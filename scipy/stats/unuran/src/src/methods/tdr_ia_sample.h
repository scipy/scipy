/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_ia_sample.c                                              *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection                                *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Sampling routines for variant with proportional squeezes and         *
 *      immediate acceptance.                                                *
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
_unur_tdr_ia_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (immediate acceptance)                         */
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
  int use_ia;
  double U, V, X;
  double fx, hx, Thx;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);

  if (GEN->iv == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"empty generator object");
    return INFINITY;
  } 

  /* main URNG */
  urng = gen->urng;

  while (1) {

    /* sample from U(0,1) */
    U = _unur_call_urng(urng);

    /* look up in guide table and search for segment */
    iv =  GEN->guide[(int) (U * GEN->guide_size)];
    U *= GEN->Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum;    /* result: U in (-A_hat,0) */

    /* check for region of immediate acceptance */
    if (U >= - iv->sq * iv->Ahat) {
      /* region of immediate acceptance */
      U /= iv->sq;
      use_ia = 1;
    }
    else {
      /* rejection from region between hat and squeeze */
      U = (U + iv->sq * iv->Ahat) / (1. - iv->sq);
      use_ia = 0;
    }
    /* result: U in (-A_hat,0) */

    /* U in (-A_hatl, A_hatr) */
    U += iv->Ahatr;

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
	X = iv->x + U /iv->fx;
      else {
	U *= iv->Tfx; /* avoid one multiplication */
	X = iv->x + (iv->Tfx * U) / (1. - iv->dTfx * U);  
	/* It cannot happen, that the denominator becomes 0 ! */
      }
      break;

    case TDR_VAR_T_POW:
      /** TODO **/
      return 1.;
      break;

    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 1.;

    } /* end switch */

    /* immedate acceptance */
    if (use_ia)
      return X;

    /* evaluate hat at X */
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      hx = iv->fx * exp(iv->dTfx*(X - iv->x)); break;
    case TDR_VAR_T_SQRT:
      Thx = iv->Tfx + iv->dTfx * (X - iv->x);     /* transformed hat at X */ 
      hx = 1./(Thx*Thx); break;
    case TDR_VAR_T_POW:
    default:
      /** TODO **/
      return 1.;
    } /* end switch */

    /* from now on we use the auxilliary generator
       (it can be the same as the main generator) */
    urng = gen->urng_aux;

    /* rejection from region between hat and (proportional) squeeze */
    V = _unur_call_urng(urng);

    /* get uniform random number between squeeze(X) and hat(X) */
    V = (iv->sq + (1 - iv->sq) * V) * hx;

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
  }

} /* end of _unur_tdr_ia_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_tdr_ia_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify results (immediate acceptance)      */
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
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng;             /* pointer to uniform random number generator */
  struct unur_tdr_interval *iv;
  int use_ia;
  double U, V, X;
  double fx, hx, Thx, sqx;
  int error = 0;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);

  if (GEN->iv == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"empty generator object");
    return INFINITY;
  } 

  /* main URNG */
  urng = gen->urng;

  while (1) {

    /* sample from U(0,1) */
    U = _unur_call_urng(urng);

    /* look up in guide table and search for segment */
    iv =  GEN->guide[(int) (U * GEN->guide_size)];
    U *= GEN->Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum;    /* result: U in (-A_hat,0) */

    /* check for region of immediate acceptance */
    if (U >= - iv->sq * iv->Ahat) {
      /* region of immediate acceptance */
      U /= iv->sq;
      use_ia = 1;
    }
    else {
      /* rejection from region between hat and squeeze */
      U = (U + iv->sq * iv->Ahat) / (1. - iv->sq);
      use_ia = 0;
    }
    /* result: U in (-A_hat,0) */

    /* U in (-A_hatl, A_hatr) */
    U += iv->Ahatr;

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
	X = iv->x + U /iv->fx;
      else {
	U *= iv->Tfx; /* avoid one multiplication */
	X = iv->x + (iv->Tfx * U) / (1. - iv->dTfx * U);  
	/* It cannot happen, that the denominator becomes 0 ! */
      }
      break;

    case TDR_VAR_T_POW:
      /** TODO **/
      return 1.;
      break;

    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 1.;

    } /* end switch */

    /* evaluate hat at X */
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      hx = iv->fx * exp(iv->dTfx*(X - iv->x)); break;
    case TDR_VAR_T_SQRT:
      Thx = iv->Tfx + iv->dTfx * (X - iv->x);     /* transformed hat at X */ 
      hx = 1./(Thx*Thx); break;
    case TDR_VAR_T_POW:
    default:
      /** TODO **/
      return 1.;
    } /* end switch */

    /* evaluate PDF at X */
    fx = PDF(X);

    /* evaluate squeeze */
    sqx = iv->sq*hx;

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

    /* immedate acceptance */
    if (use_ia)
      return X;

    /* from now on we use the auxilliary generator
       (it can be the same as the main generator) */
    urng = gen->urng_aux;

    /* rejection from region between hat and (proportional) squeeze */
    V = _unur_call_urng(urng);

    /* get uniform random number between squeeze(X) and hat(X) */
    V = (iv->sq + (1 - iv->sq) * V) * hx;

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
  }

} /* end of _unur_tdr_ia_sample_check() */

/*---------------------------------------------------------------------------*/
