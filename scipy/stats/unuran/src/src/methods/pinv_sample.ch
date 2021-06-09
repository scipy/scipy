/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      pinv_sample.c                                                *
 *                                                                           *
 *   Sampling routines.                                                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008 Wolfgang Hoermann and Josef Leydold                  *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/**  Sampling routines                                                      **/
/*****************************************************************************/

double
_unur_pinv_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
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
  double U,X;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_PINV_GEN,INFINITY);

  /* sample from U(0,1) */
  U = _unur_call_urng(gen->urng);

  /* compute inverse CDF */
  X = _unur_pinv_eval_approxinvcdf(gen,U);

  if (X<DISTR.trunc[0]) return DISTR.trunc[0];
  if (X>DISTR.trunc[1]) return DISTR.trunc[1];

  return X;

} /* end of _unur_pinv_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_eval_approxinvcdf( const struct unur_gen *gen, double u )
     /*----------------------------------------------------------------------*/
     /* evaluate polynomial interpolation of inverse CDF at u                */
     /* (internal call)                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   u   ... argument for inverse CDF (0<=u<=1, no validation!)         */
     /*                                                                      */
     /* return:                                                              */
     /*   double (approximate inverse CDF)                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  int i;
  double x,un;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_PINV_GEN,INFINITY);

  /* rescale for range (0, Umax) */
  un = u * GEN->Umax;

  /* look up in guide table and search for interval */
  i = GEN->guide[(int)(u * GEN->guide_size)];
  while (GEN->iv[i+1].cdfi < un)
    i++;

  /* rescale for range (0, CDF(right)-CDF(left) for interval */
  un -= GEN->iv[i].cdfi;

  /* evaluate polynomial */
  x = _unur_pinv_newton_eval(un, GEN->iv[i].ui, GEN->iv[i].zi, GEN->order);

  /* return point (add left boundary point to x) */
  return (GEN->iv)[i].xi + x;

} /* end of _unur_pinv_eval_approxinvcdf() */

/*---------------------------------------------------------------------------*/

double
unur_pinv_eval_approxinvcdf( const struct unur_gen *gen, double u )
     /*----------------------------------------------------------------------*/
     /* evaluate polynomial interpolation of inverse CDF at u                */
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
  if ( gen->method != UNUR_METH_PINV ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INFINITY;
  }
  COOKIE_CHECK(gen,CK_PINV_GEN,INFINITY);

  if ( ! (u>0. && u<1.)) {
    if ( ! (u>=0. && u<=1.)) {
      _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"U not in [0,1]");
    }
    if (u<=0.) return DISTR.domain[0];
    if (u>=1.) return DISTR.domain[1];
    return u;  /* = NaN */
  }
  
  /* compute inverse CDF */
  x = _unur_pinv_eval_approxinvcdf(gen,u);

  /* validate range */
  if (x<DISTR.domain[0]) x = DISTR.domain[0];
  if (x>DISTR.domain[1]) x = DISTR.domain[1];

  return x;

} /* end of unur_pinv_eval_approxinvcdf() */

/*---------------------------------------------------------------------------*/

double
unur_pinv_eval_approxcdf( const struct unur_gen *gen, double x )
     /*----------------------------------------------------------------------*/
     /* evaluate (approximate) CDF at x.                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... argument for CDF                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   double (approximate CDF)                                           */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  if ( gen->method != UNUR_METH_PINV ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INFINITY;
  }
  COOKIE_CHECK(gen,CK_PINV_GEN,INFINITY);

  /* we need the table of CDF values */
  if ( (gen->variant & PINV_VARIANT_PDF) && GEN->aCDF == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GENERIC,"'keepcdf' not set");
    return INFINITY;
  }

  /* argument inside domain ? */
  if (x <= DISTR.domain[0]) return 0.;
  if (x >= DISTR.domain[1]) return 1.;

  /* compute CDF */
  if (gen->variant & PINV_VARIANT_PDF) {
    /* case: PDF given */
    return _unur_lobatto_eval_CDF(GEN->aCDF,x);
  }
  else {
    /* case: CDF given */
    return (CDF(x));
  }

} /* end of unur_pinv_eval_approxcdf() */

/*****************************************************************************/

int
unur_pinv_estimate_error( const UNUR_GEN *gen, int samplesize, double *max_error, double *MAE )
     /*----------------------------------------------------------------------*/
     /* Estimate maximal u-error and mean absolute error (MAE) by means of   */
     /* Monte-Carlo simulation.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   samplesize ... sample size for Monte Carlo simulation              */
     /*   max_error  ... pointer to double for storing maximal u-error       */
     /*   MAE        ... pointer to double for storing MA u-error            */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{ 
  /* check arguments */
  _unur_check_NULL(GENTYPE, gen, UNUR_ERR_NULL);  
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE);

  /* run test */
  unur_test_u_error(gen, max_error, MAE, 1.e-20, samplesize, 
		     FALSE, FALSE, FALSE, NULL);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_pinv_estimate_error() */

/*****************************************************************************/
