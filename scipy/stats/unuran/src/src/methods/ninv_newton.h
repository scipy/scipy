/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ninv_newton.ch                                               *
 *                                                                           *
 *   Routines for Newton's root finding algorithm.                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   REFERENCES:                                                             *
 *   [1] Neumaier A. (to be published): Introduction to numerical analysis,  *
 *       Cambridge University Press                                          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  For a given U find X with F(X) - U = 0.                                  *
 *                                                                           *
 *  Method: damped Newton method.                                            *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/**  Newton method                                                          **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/* maximal number of steps to leave flat region */
#define MAX_FLAT_COUNT  (40)

/*---------------------------------------------------------------------------*/

double
_unur_ninv_newton( const struct unur_gen *gen, double U )
     /*----------------------------------------------------------------------*/
     /* sample from generator (use Newton's method)                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*     U   ... random number (uniform distribution)                     */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double x;           /* point for netwon-iteration                   */
  double fx;          /* cdf at x                                     */
  double dfx;         /* pdf at x                                     */
  double fxabs;       /* absolute valuo of fx                         */
  double xtmp, fxtmp; /* temprary variables for x and fx              */
  double xold;        /* remember last values for stopping criterion  */
  double fxtmpabs;    /* fabs of fxtmp                                */
  double damp;        /* damping factor                               */
  double step;        /* helps to escape from flat regions of the cdf */
  int i;              /* counter for for-loop, index                  */
  int flat_count;     /* counter of steps in flat region              */
  double rel_u_resolution; /* relative u resolution                   */
  int x_goal, u_goal; /* whether precision goal is reached            */

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_NINV_GEN,INFINITY);

  /* compute relative u resolution */
  rel_u_resolution = ( (GEN->u_resolution > 0.) ? 
                       (GEN->Umax - GEN->Umin) * GEN->u_resolution :
                       INFINITY );

  /* -- 1. initialize starting interval -- */

  if (GEN->table_on) {
    /* -- 1a. use table -- */

    /* 0 <= i <= table_size-2  */
    if ( _unur_FP_same(GEN->CDFmin,GEN->CDFmax) ) {
      /* CDF values in table too close, so we use median point since */
      /* there is no difference between CDF values.                  */
      i = GEN->table_size/2;
    }
    else {
      i = (int) ( GEN->table_size * (U - GEN->CDFmin) / (GEN->CDFmax - GEN->CDFmin) );
      if (i<0) i = 0;
      else if (i > GEN->table_size - 2) i = GEN->table_size - 2;
    }

    if (_unur_FP_is_infinity(GEN->table[i+1])) {
      x  = GEN->table[i];
      fx = GEN->f_table[i];
    }
    else {
      x  = GEN->table[i+1];
      fx = GEN->f_table[i+1];
    }

  }

  else { 
    /* 1b. -- no table available -- */
    x  = GEN->s[0];
    fx = GEN->CDFs[0];
  }

  /* -- 1c. check for boundary of truncated domain -- */

  if ( x < DISTR.trunc[0] ){
    x  = DISTR.trunc[0];
    fx = GEN->Umin;    /* = CDF(x) */
  }
  else if ( x > DISTR.trunc[1] ){
    x  = DISTR.trunc[1];
    fx = GEN->Umax;    /* = CDF(x) */
  }

  /* -- 1z. compute values for starting point -- */

  fx   -= U;
  dfx   = PDF(x);
  fxabs = fabs(fx);
  xold  = x;    /* there is no old value yet */

  damp = 2.;        /* to be halved at least once */  
  step = 1.;

  /* -- 2. Newton iteration -- */

  for (i=0; i < GEN->max_iter; i++) {

    flat_count = 0;
    while (_unur_iszero(dfx)) {   /* function flat at x */
      /* printf("step: %g, x: %g, fx: %g, dfx: %g\n",step, x, fx, dfx); */

      if (_unur_iszero(fx))  /* exact hit -> leave while-loop */
	break; 

      if (fx > 0.) {         /* try another x */
        xtmp = x - step; 
	xtmp = _unur_max( xtmp, DISTR.domain[0] );
      }
      else {
        xtmp  = x + step;
	xtmp = _unur_min( xtmp, DISTR.domain[1] );
      }

      fxtmp    = CDF(xtmp) - U;
      fxtmpabs = fabs(fxtmp);

      if ( fxtmpabs < fxabs ) {       /* improvement, update x               */
	/* printf("fxabs: %g tmpabs: %g\n", fxabs, fxtmpabs); */
        step = 1.;     /* set back stepsize */
        x     = xtmp;
        fx    = fxtmp;
      }
      else if ( fxtmp*fx < 0. ) {     /* step was too large, don't update x  */
        step /= 2.;                      
      } 
      else {                          /* step was too short, update x        */
        step *= 2.;    
        x     = xtmp;
        fx    = fxtmp;
      }  

      dfx   = PDF(x);
      fxabs = fabs(fx);

      if (flat_count < MAX_FLAT_COUNT)
	flat_count++;
      else {
	_unur_error(gen->genid,UNUR_ERR_GEN_SAMPLING,
		    "Newton's method cannot leave flat region");
	x = _unur_max( x, DISTR.trunc[0]);
	x = _unur_min( x, DISTR.trunc[1]);
	return x;
      }
    }   /* end of while-loop, (leaving flat region) */
    
    step = 1.;   /* set back stepsize */

    if (_unur_iszero(fx))  /* exact hit -> finished */
      break;


    if (_unur_isfinite(dfx)) {
      do {    /* newton-step  (damped if nececcary) */
	damp /= 2.;
	xtmp = x - damp * fx/dfx;
	/* make sure that new point is inside (truncated) domain */
	xtmp = _unur_min( xtmp, DISTR.trunc[1] );
	xtmp = _unur_max( xtmp, DISTR.trunc[0] );
	fxtmp = CDF(xtmp) - U;
      } while (fabs(fxtmp) > fxabs * (1.+UNUR_SQRT_DBL_EPSILON));   /* no improvement */
    }
    else {
      /* we cannot use Newton's rule if the derivative is not finite. */
      /* this happens when we hit a pole of the PDF.                  */
      /* use a bisection step instead.                                */
      xtmp = 0.5*(x + xold);
      fxtmp = CDF(xtmp) - U;
    }
    
    /* updation variables according to newton-step      */
    damp  = 2.;       /* set back factor for damping    */
    xold  = x;        /* remember last value of x       */
    x     = xtmp;     /* update x                       */
    fx    = fxtmp;    /* update function value at x     */
    dfx   = PDF(x);   /* update derivative sof fx at x  */
    fxabs = fabs(fx); /* update absolute value of fx    */
 
    /* -- 2z. check stopping criterions -- */

    if ( GEN->x_resolution > 0. ) {
      /* check x-error */
      /* we use a combination of absolute and relative x-error: */
      /*    x-error < x-resolution * fabs(x) + x-resolution^2   */
      if ( _unur_iszero(fx) ||                            /* exact hit */ 
           fabs(x-xold) < GEN->x_resolution * (fabs(x) + GEN->x_resolution) ) {
 	x_goal = TRUE;
      }
      else
        x_goal = FALSE;
    }
    else {
      /* no check */
      x_goal = TRUE;
    }

    if ( GEN->u_resolution > 0. ) {
      /* check u-error */
      /* (we use a slightly smaller maximal tolerated error than given by user) */
      if ( fabs(fx) < 0.9 * rel_u_resolution ) {    /* relative u resolution */
      	u_goal = TRUE;
      }
      else if ( _unur_FP_same(xold, x) ) {
        /* sharp peak or pole */
        _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
                      "sharp peak or pole: accuracy goal in u cannot be reached");
        u_goal = TRUE;
      }
      else
        u_goal = FALSE;

    }
    else {
      u_goal = TRUE;
    }

    /* goal reached ? */
    if (x_goal && u_goal)
      /*finished*/
      break;
  }  /* end of for-loop  (MAXITER reached -> finished) */


  if (i >= GEN->max_iter)
    _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
		  "max number of iterations exceeded: accuracy goal might not be reached");

  /* make sure that result is within boundaries of (truncated) domain */
  x = _unur_max( x, DISTR.trunc[0]);
  x = _unur_min( x, DISTR.trunc[1]);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file (in case error) */
  if (gen->debug & NINV_DEBUG_SAMPLE)
    _unur_ninv_debug_sample(gen, U, x, fx, i);
#endif
  
  return x;

} /* end of _unur_ninv_sample_newton() */

/*---------------------------------------------------------------------------*/

#undef MAX_FLAT_COUNT

/*---------------------------------------------------------------------------*/
