/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ninv_sample.ch                                               *
 *                                                                           *
 *   Sampling routines.                                                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/**  Sampling routines                                                      **/
/*****************************************************************************/

double 
_unur_ninv_sample_newton( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (use newtons method)                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*----------------------------------------------------------------------*/
{
  return _unur_ninv_newton( gen, 
         GEN->Umin + (_unur_call_urng(gen->urng)) * (GEN->Umax - GEN->Umin) );
}

/*---------------------------------------------------------------------------*/

double
_unur_ninv_sample_regula( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (use regula falsi)                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*----------------------------------------------------------------------*/
{
  return _unur_ninv_regula( gen, 
         GEN->Umin + (_unur_call_urng(gen->urng)) * (GEN->Umax - GEN->Umin) );
} /* end of _unur_ninv_sample_regula() */

/*---------------------------------------------------------------------------*/

double
_unur_ninv_sample_bisect( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (use bisection method)                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*----------------------------------------------------------------------*/
{
  return _unur_ninv_bisect( gen, 
         GEN->Umin + (_unur_call_urng(gen->urng)) * (GEN->Umax - GEN->Umin) );
} /* end of _unur_ninv_sample_bisect() */

/*---------------------------------------------------------------------------*/

double
unur_ninv_eval_approxinvcdf( const struct unur_gen *gen, double u )
     /*----------------------------------------------------------------------*/
     /* get approximate value of inverse CDF at u approximately              */
     /* (user call)                                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   u   ... argument for inverse CDF (0<=u<=1)                         */
     /*                                                                      */
     /* return:                                                              */
     /*   double (approximate inverse CDF)                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double x;

  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  if ( gen->method != UNUR_METH_NINV ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INFINITY; 
  }
  COOKIE_CHECK(gen,CK_NINV_GEN,INFINITY);

  if ( ! (u>0. && u<1.)) {
    if ( ! (u>=0. && u<=1.)) {
      _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"U not in [0,1]");
    }
    if (u<=0.) return DISTR.domain[0];
    if (u>=1.) return DISTR.domain[1];
    return u;  /* = NaN */
  }
  
  /* compute inverse CDF */
  switch (gen->variant) {
  case NINV_VARFLAG_NEWTON:
    x = _unur_ninv_newton(gen,u);
    break;
  case NINV_VARFLAG_BISECT:
    x = _unur_ninv_bisect(gen,u);
    break;
  case NINV_VARFLAG_REGULA:
  default:
    x = _unur_ninv_regula(gen,u);
    break;
  }

  /* validate range */
  if (x<DISTR.domain[0]) x = DISTR.domain[0];
  if (x>DISTR.domain[1]) x = DISTR.domain[1];

  return x;

} /* end of unur_hinv_eval_approxinvcdf() */

/*---------------------------------------------------------------------------*/
